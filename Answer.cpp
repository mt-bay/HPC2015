//------------------------------------------------------------------------------
/// @file
/// @brief    HPCAnswer.hpp の実装 (解答記述用ファイル)
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2015 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください

//------------------------------------------------------------------------------

#include "HPCAnswer.hpp"
#include "HPCMath.hpp"

#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cstdint>

/// プロコン問題環境を表します。
namespace hpc {
    class Path;

    void choose_item(const Field& _field, const ItemCollection& _items);
    std::vector<Action> gen_route(const Field& _field, const ItemCollection& _items, const std::list<int>& _load_index);
    std::list<int> convex_hull(const Field& _field, const ItemCollection& _items, const std::vector<int>& _load_index);
    std::vector<int> sort_position_y(const Field& _field, const ItemCollection& _items, const std::vector<int>& _target);
    const Pos get_pos(const Field& _field, const ItemCollection& _items, const int _index);
    int get_weight(const ItemCollection& _items, const int _index);
    int distance(const Pos& _a, const Pos& _b);
    int get_delivery_cost(const std::list<Path>& _path, const int _weight, std::list<Path>::iterator _branch_timing, const Path& _add_path_start_to_cp, const Path& _add_path_cp_to_goal);
    int get_delivery_cost(const int _weight, const Path& _add_path_start_to_cp, const Path& _add_path_cp_to_goal);

    //------------------------------------------------------------------------------
    /// @brief startからendまでの経路を示します
    class Path
    {
    public:
        int start;
        int goal;

        int load_weight;
        std::list<Action> route;

    private:
        class AStar_dat
        {
        private:
            inline int get_huries(const Pos& _goal, const Pos& _here)
                const
            {
                return Math::Abs(_goal.x - _here.x) + Math::Abs(_goal.y - _here.y);
            }

        public:
            enum State
            {
                NONE,
                OPEN,
                CLOSE
            };

            AStar_dat*    m_parent;
            Pos   m_pos;

            State m_state;
            int   m_cost;
            int   m_huries;

            // デフォルトコンストラクタ
            AStar_dat()
                :
                m_parent(nullptr),
                m_state (NONE)   ,
                m_cost  (0)      ,
                m_huries(0)
            {
            }

            // コピーコンストラクタ
            AStar_dat(const AStar_dat& _orig)
                :
                m_parent(_orig.m_parent),
                m_state (_orig.m_state) ,
                m_cost  (_orig.m_cost)  ,
                m_huries(_orig.m_huries)
            {
            }

            // ノードを開く
            // @return ノードがOPENになったか
            bool open(const Pos& _here, const Pos& _goal, AStar_dat* _parent = nullptr, const bool _is_wall = false)
            {
                if(_is_wall)
                {
                    m_state = CLOSE;
                    return false;
                }

                m_state  = OPEN;
                m_parent = _parent;
                m_pos    = _here;
                m_cost   = (m_parent == nullptr ? 0 : m_parent->m_cost) + 1;
                m_huries  = get_huries(_goal, _here);

                return true;
            }

            inline int get_score()
                const
            {
                return m_cost + m_huries;
            }
        };

        // ゴールからスタートまでのノードがOPENになっているデータからPATH導出
        std::list<Action> gen_path(std::vector< std::vector<AStar_dat> >& _data, const Pos& _goal)
            const
        {
            std::list<Action> buf_return;

            AStar_dat* buf_data = &_data.at(_goal.x).at(_goal.y), *buf_parent = buf_data->m_parent;
            while(buf_parent != nullptr)
            {
                if     (buf_data->m_pos.x - buf_parent->m_pos.x ==  1)
                {
                    buf_return.push_front(Action_MoveRight);
                }
                else if(buf_data->m_pos.x - buf_parent->m_pos.x == -1)
                {
                    buf_return.push_front(Action_MoveLeft);
                }
                else if(buf_data->m_pos.y - buf_parent->m_pos.y ==  1)
                {
                    buf_return.push_front(Action_MoveUp);
                }
                else if(buf_data->m_pos.y - buf_parent->m_pos.y == -1)
                {
                    buf_return.push_front(Action_MoveDown);
                }

                //std::cout << "\ttarget : (" << buf_data->m_pos.x << ", " << buf_data->m_pos.y << "), parent : (" << buf_parent->m_pos.x << ", " << buf_parent->m_pos.y << "), act : " << buf_return.front() << std::endl;
                buf_data   = &_data.at(buf_parent->m_pos.x).at(buf_parent->m_pos.y);
                buf_parent = buf_parent->m_parent;
            }

            return buf_return;
        }

        // A*探索
        std::list<Action> a_star(const Field& _field, const Pos& _start, const Pos& _goal)
            const
        {
            if(_goal == _start)
            {
                return std::list<Action>();
            }

            //std::cout << "AStar((" << _start.x << ", " << _start.y << ")->(" << _goal.x << ", " << _goal.y <<")" << std::endl;

            std::vector< std::vector<AStar_dat> > data(_field.width(), std::vector<AStar_dat>(_field.height(), AStar_dat()));

            std::vector<AStar_dat*> open_nodes;  //Open状態のノードリスト
            std::list<AStar_dat*>   close_nodes; //Close状態のノードリスト

            //初期ノード生成
            data.at(_start.x).at(_start.y).open(_start, _goal);
            open_nodes.push_back(&data.at(_start.x).at(_start.y));

            Pos    buf_pos, buf_move;

            AStar_dat* ref_buf_data;
            //探索
            int open_num = 0;
            while(true)
            {
                // スコア最小のノードを導出
                int buf_index = 0, buf_score = open_nodes.front()->get_score();
                for(size_t s = 1, cond_s = open_nodes.size(); s < cond_s; ++s)
                {
                    if(buf_score > open_nodes.at(s)->get_score())
                    {
                        buf_index = s;
                        buf_score = open_nodes.at(s)->get_score();
                    }
                }
                buf_pos = open_nodes.at(buf_index)->m_pos;
                // スコア最小のノードの上下左右を探索する(OPENにする)
                for(int i = 0; i < 4; ++i)
                {
                    buf_move = buf_pos.move((Action)i);
                    ref_buf_data = &data.at(buf_move.x).at(buf_move.y);
                    if(ref_buf_data->m_state == AStar_dat::NONE)
                    {
                        if(ref_buf_data->open(buf_move, _goal, open_nodes.at(buf_index), _field.isWall(buf_move)))
                        {
                            open_nodes.push_back(ref_buf_data);
                        }
                        ++open_num;
                        if(ref_buf_data->m_pos == _goal)
                        {
                            return gen_path(data, _goal);
                        }
                    }
                }
                open_nodes.at(buf_index)->m_state = AStar_dat::CLOSE;
                open_nodes.erase(std::next(open_nodes.begin(), buf_index));
            }
        }

    public:
        // A*を用いてstartからgoalまでの経路導出
        Path(const Field& _field ,const ItemCollection& _item ,const int _start, const int _goal)
            :
            start      (_start),
            goal       (_goal),
            load_weight(get_weight(_item, _goal)),
            route      (a_star(_field, get_pos(_field, _item, _start), get_pos(_field, _item, _goal)))
        {
        }

        Path()
            :
            start      (0),
            goal       (0),
            load_weight(0),
            route      (std::list<Action>())
        {
        }

        //コピーコンストラクタ
        Path(const Path& _orig)
            :
            start      (_orig.start),
            goal       (_orig.goal) ,
            load_weight(_orig.load_weight),
            route      (_orig.route)
        {

        }

        Path& operator=(const Path& _rhs)
        {
            start       = _rhs.start;
            goal        = _rhs.goal;
            load_weight = _rhs.load_weight;
            route       = _rhs.route;

            return *this;
        }

        inline size_t cost()
            const
        {
            return route.size();
        }

        inline std::list<Action> get_route()
            const
        {
            return route;
        }

        inline bool start_is_office(const int item_count)
            const
        {
            return start < 0 || item_count <= start;
        }
    };

    namespace 
    {
        std::vector<Action> action[Parameter::PeriodCount];
        std::list<int>      load_index[Parameter::PeriodCount];
        int    period;
        size_t act_index;
    };


    //------------------------------------------------------------------------------
    /// 各ステージ開始時に呼び出されます。
    ///
    /// ここで、各ステージに対して初期処理を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    void Answer::Init(const Stage& aStage)
    {
        choose_item(aStage.field(), aStage.items());
        //std::cout << "num. of item : " << aStage.items().count() << std::endl;
        
    }

    //------------------------------------------------------------------------------
    /// 各配達時間帯開始時に呼び出されます。
    ///
    /// ここで、この時間帯に配達する荷物をトラックに積み込みます。
    /// どの荷物をトラックに積み込むかを決めて、引数の aItemGroup に対して
    /// setItem で荷物番号を指定して積み込んでください。
    ///
    /// @param[in] aStage 現在のステージ。
    /// @param[in] aItemGroup 荷物グループ。
    void Answer::InitPeriod(const Stage& aStage, ItemGroup& aItemGroup)
    {
        //std::cout << "\tperiod : " << period << ", item : ";
        for(std::list<int>::const_iterator it = load_index[aStage.period()].cbegin(), cond_it = load_index[aStage.period()].cend(); it != cond_it; ++it)
        {
            //std::cout << "(" << *it << ", " << aStage.items()[*it].period() << "), ";
            aItemGroup.addItem(*it);
        }
        
        period    = aStage.period();
        act_index = 0;
        //std::cout << "size of action : " << action[period].size() << std::endl;
        return;
    }

    //------------------------------------------------------------------------------
    /// 各ターンでの動作を返します。
    ///
    /// @param[in] aStage 現在ステージの情報。
    ///
    /// @return これから行う動作を表す Action クラス。
    Action Answer::GetNextAction(const Stage& aStage)
    {
        Action next_act = action[period].at(act_index);

        //std::cout << "\t\tturn : " << act_index << ", truck : (" << aStage.truck().pos().x << ", " << aStage.truck().pos().y << "), next act : " << next_act << std::endl;
        ++act_index;
        return next_act;
    }

    //------------------------------------------------------------------------------
    /// 各配達時間帯終了時に呼び出されます。
    ///
    /// ここで、この時間帯の終了処理を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    /// @param[in] aStageState 結果。Playingならこの時間帯の配達完了で、それ以外なら、何らかのエラーが発生した。
    /// @param[in] aCost この時間帯に消費した燃料。エラーなら0。
    void Answer::FinalizePeriod(const Stage& aStage, StageState aStageState, int aCost)
    {
        //std::cout << std::endl;
    }

    //------------------------------------------------------------------------------
    /// 各ステージ終了時に呼び出されます。
    ///
    /// ここで、各ステージに対して終了処理を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    /// @param[in] aStageState 結果。Completeなら配達完了で、それ以外なら、何らかのエラーが発生した。
    /// @param[in] aScore このステージで獲得したスコア。エラーなら0。
    void Answer::Finalize(const Stage& aStage, StageState aStageState, int aScore)
    {
        for(size_t s = 0; s < Parameter::PeriodCount; ++s)
        {
            load_index[s].clear();
            action[s].clear();
        }
    }


    //------------------------------------------------------------------------------
    //以下，追加関数
    
    //積載物選択
    //4ピリオドの荷物全部決定する
    void choose_item(const Field& _field, const ItemCollection& _items)
    {
        static int      capacity = Parameter::TruckWeight + Parameter::TruckWeightCapacity;
        std::list<int>  free_item[Parameter::ItemWeightMax];
        std::list<Path> delivery_cost[Parameter::PeriodCount];
        int             inital_weight[Parameter::PeriodCount] = {Parameter::TruckWeight, Parameter::TruckWeight, Parameter::TruckWeight, Parameter::TruckWeight};
        int             buf_weight;

        int min_diff = std::numeric_limits<int32_t>::max();
        std::list<Path>::iterator min_it;
        Path min_start_to_s, min_s_to_goal;
        Path buf_start_to_s, buf_s_to_goal;
        int buf_cost = get_delivery_cost(delivery_cost[period], inital_weight[period], delivery_cost[period].end(), buf_start_to_s, buf_s_to_goal);
        int buf_diff;
        int period;

        //各ピリオドで配達することを指定されているものはとりあえず載せる
        for(int s = 0, cond_s = _items.count(); s < cond_s; ++s)
        {
            if(_items[s].period() == -1)
            {
                free_item[_items[s].weight() - 1].push_back(s);
            }
            else
            {
                min_diff = std::numeric_limits<int32_t>::max();
                period = _items[s].period();

                if(load_index[period].empty())
                {
                    delivery_cost[period].push_back(Path(_field, _items, _items.count(), s));
                    inital_weight[period] += delivery_cost[period].back().load_weight;
                    delivery_cost[period].push_back(Path(_field, _items, s, _items.count()));
                    inital_weight[period] += delivery_cost[period].back().load_weight;
                }
                else
                {
                    buf_cost = get_delivery_cost(delivery_cost[period], inital_weight[period], delivery_cost[period].end(), buf_start_to_s, buf_s_to_goal);
                    buf_weight = inital_weight[period] + get_weight(_items, s);

                    for(std::list<Path>::iterator it = delivery_cost[period].begin(), cond_it = delivery_cost[period].end(); it != cond_it; ++it)
                    {
                        if(min_diff > (buf_diff = get_delivery_cost(delivery_cost[period], buf_weight, it, (buf_start_to_s = Path(_field, _items, it->start, s)), (buf_s_to_goal = Path(_field, _items, s, it->goal))) - buf_cost))
                        {
                            min_diff = buf_diff;
                            min_it = it;
                            min_start_to_s = buf_start_to_s;
                            min_s_to_goal  = buf_s_to_goal;
                        }
                    }
                    min_it = delivery_cost[period].erase(min_it);
                    min_it = delivery_cost[period].insert(min_it, min_s_to_goal);
                    min_it = delivery_cost[period].insert(min_it, min_start_to_s);
                    inital_weight[period] = buf_weight;
                }
                load_index[period].push_back(s);

            }
        }
        
        //時間帯指定なしの荷物を，重いものから載せていく
        //載せるピリオドが決まっているものから最近傍点を導出し，その点と同じピリオドに載せることにする
        int min_period, target_weight = Parameter::ItemWeightMax;
        bool empty_is_searched;
        for(int i = Parameter::ItemWeightMax - 1; i >= 0; --i)
        {
            while(!free_item[i].empty())
            {
                empty_is_searched = false;
                min_period = 0;
                min_diff   = std::numeric_limits<int32_t>::max();
                min_it     = delivery_cost[min_period].end();

                for(int t = 0; t < Parameter::PeriodCount; ++t)
                {
                    if(capacity < (buf_weight = inital_weight[t] + target_weight))
                    {
                        continue;
                    }

                    if(load_index[t].empty() && !empty_is_searched)
                    {
                        if(min_diff > (buf_diff = get_delivery_cost(buf_weight, (buf_start_to_s = Path(_field, _items, _items.count(), free_item[i].front())), (buf_s_to_goal = Path(_field, _items, free_item[i].front(), _items.count())))))
                        {
                            min_diff       = buf_diff;
                            min_period     = t;
                            min_it         = delivery_cost[t].end();
                            min_start_to_s = buf_start_to_s;
                            min_s_to_goal  = buf_s_to_goal;
                        }
                        empty_is_searched = true;
                    }
                    else
                    {
                        buf_cost = get_delivery_cost(delivery_cost[t], inital_weight[t], delivery_cost[t].end(), buf_start_to_s, buf_s_to_goal);
                        for(std::list<Path>::iterator it = delivery_cost[t].begin(), cond_it = delivery_cost[t].end(); it != cond_it; ++it)
                        {
                            if(min_diff > (buf_diff = get_delivery_cost(delivery_cost[t], buf_weight, it, (buf_start_to_s = Path(_field, _items, it->start, free_item[i].front())), (buf_s_to_goal = Path(_field, _items, free_item[i].front(), it->goal ))) - buf_cost))
                            {
                                min_diff       = buf_diff;
                                min_period     = t;
                                min_it         = it;
                                min_start_to_s = buf_start_to_s;
                                min_s_to_goal  = buf_s_to_goal;
                            }
                        }
                    }
                }
                
                load_index[min_period].push_back(free_item[i].front());
                if(min_it != delivery_cost[min_period].end())
                {
                    min_it = delivery_cost[min_period].erase(min_it);
                }
                min_it = delivery_cost[min_period].insert(min_it, min_s_to_goal);
                min_it = delivery_cost[min_period].insert(min_it, min_start_to_s);
                inital_weight[min_period] += target_weight;

                free_item[i].pop_front();
            }
            --target_weight;
        }


        //巡回ルートのアクションを各ピリオド用のグローバル変数に投げる
        for(int i = 0; i < Parameter::PeriodCount; ++i)
        {
            for(std::list<Path>::iterator it = delivery_cost[i].begin(), cond_it = delivery_cost[i].end(); it != cond_it; ++it)
            {
                //std::cout << "add in period " << i << "(size : " << it->cost() << "), (" << get_pos(_field, _items, it->start).x << ", " << get_pos(_field, _items, it->start).y << ")->(" << get_pos(_field, _items, it->goal).x << ", " << get_pos(_field, _items, it->goal).y << ")" << std::endl;
                action[i].insert(action[i].end(), it->route.begin(), it->route.end());
            }
        }

        return;
    }


    //ピリオド単位の輸送コスト導出
    //_add_path_index = end()ならば追加なしの経路コスト導出, それ以外ならばしかるべき場所で_add_pathを経由した経路のコスト導出
    int get_delivery_cost(const std::list<Path>& _path, const int _weight, std::list<Path>::iterator _branch_timing, const Path& _add_path_start_to_cp, const Path& _add_path_cp_to_goal)
    {
        int buf_return = 0, buf_weight = _weight;

        for(std::list<Path>::const_iterator it = _path.cbegin(), cond_it = _path.cend(); it != cond_it; ++it)
        {
            if(it == _branch_timing)
            {
                buf_return += _add_path_start_to_cp.cost() * _weight;
                buf_weight -= _add_path_start_to_cp.load_weight;

                buf_return += _add_path_cp_to_goal.cost() * _weight;
                buf_weight -= _add_path_cp_to_goal.load_weight;
            }
            else
            {
                buf_return += it->cost() * _weight;
                buf_weight -= it->load_weight;
            }
        }
        return buf_return;
    }


    int get_delivery_cost(const int _weight, const Path& _add_path_start_to_cp, const Path& _add_path_cp_to_goal)
    {
        int buf_return = 0, buf_weight = _weight;

        buf_return += _add_path_start_to_cp.cost() * _weight;
        buf_weight -= _add_path_start_to_cp.load_weight;

        buf_return += _add_path_cp_to_goal.cost() * _weight;

        return buf_return;
    }


    const Pos get_pos(const Field& _field, const ItemCollection& _items, const int _index)
    {
        return (int)0 <= _index && _index < _items.count() ? _items[_index].destination() : _field.officePos();
    }


    int get_weight(const ItemCollection& _items, const int _index)
    {
        return (int)0 <= _index && _index < _items.count() ? _items[_index].weight() : 0;
    }


    int distance(const Pos& _a, const Pos& _b)
    {
        return Math::Abs(_b.x - _a.x) + Math::Abs(_b.y - _a.y);
    }
}

//------------------------------------------------------------------------------
// EOF
